const render = (thet,P0,h,t1,t2,i) => {
    const app = document.querySelector('#app');
    let class_name = `data-${i}`;
    const container = document.createElement('div');
    container.classList.add('app__box');
    app.insertAdjacentElement('beforeend', container);
    const block = `
    <h2 class="app__title">Расчеты №${i}</h2>
    <ul class="app__fix"></ul>
    <div class="app__table">
    <table>
    <thead>
    <tr>
    <th>t, с</th>
    <th>m, кг</th>
    <th>Vx, м/с</th>
    <th>Vy, м/с</th>
    <th>x, м</th>
    <th>y, м </th>
    <th>V, м/с </th>
    <th>r, м</th>
    <th>thet, рад</th>
    <th>THET , рад</th>
    <th>alfa, рад</th>
    <th>fi, рад</th>
    <th>THETc, рад</th>
    </tr>
    </thead>
    </table>
    <div class="app__table-body">
    <table>
    <tbody id="${class_name}">
    </tbody>
    </table>
    </div>
    </div>`
    container.insertAdjacentHTML('beforeend', block);
    const fixDiv = container.querySelector('.app__fix');
    
    fixDiv.innerHTML = `
    <li><b>P = </b>${P0} [Н]</li>
    <li><b>h = </b>${h} [м]</li>
    <li><b>t1 = </b>${t1} [c]</li>
    <li><b>t2 = </b>${t2} [c]</li>
    <li><b>thet_torch = </b>${thet[0]} [рад]</li>
    <li><b>thet_2 = </b>${thet[1]} [рад]</li>
    `
}
const html = (t, el, i) => {
    let class_name = `data-${i}`;
    const tr_iter = (t, el) => {
        return `
        <td class="app__ceil">${t.toFixed(5)}</td>
        <td class="app__ceil">${el[0].toFixed(5)}</td>
        <td class="app__ceil">${el[1].toFixed(5)}</td>
        <td class="app__ceil">${el[2].toFixed(5)}</td>
        <td class="app__ceil">${el[3].toFixed(5)}</td>
        <td class="app__ceil">${el[4].toFixed(5)}</td>
        <td class="app__ceil">${el[5].toFixed(5)}</td>
        <td class="app__ceil">${el[6].toFixed(5)}</td>
        <td class="app__ceil">${el[7].toFixed(5)}</td>
        <td class="app__ceil">${el[8].toFixed(5)}</td>
        <td class="app__ceil">${el[9].toFixed(5)}</td>
        <td class="app__ceil">${el[10].toFixed(5)}</td>
        <td class="app__ceil">${el[11].toFixed(5)}</td>
        `
    }
    
    const tbody = document.querySelector(`#${class_name}`);
    const tr = document.createElement('tr')
    tr.insertAdjacentHTML('beforeend', tr_iter(t, el))
    tbody.insertAdjacentElement('beforeend', tr)
}

const consoleLog = (t, el) => {
    console.log(`${t.toFixed(5)} ${el[0].toFixed(5)} ${el[1].toFixed(5)} ${el[2].toFixed(5)} ${el[3].toFixed(5)} ${el[4].toFixed(5)} ${el[5].toFixed(5)} ${el[6].toFixed(5)} ${el[7].toFixed(5)} ${el[8].toFixed(5)} ${el[9].toFixed(5)} ${el[10].toFixed(5)} ${el[11].toFixed(5)}`);
}

const printTable = (thet, h, t1, t2, P0) => {
    render(thet, P0, h, t1, t2, i);
    let dt = param.step;
    let t = param.t0;
    let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0, param.thet, param.THET, param.alfa, param.fi, param.THETc];
    html(t, el, i);

    let t_prev = t;
    let el_prev = el.slice();
    while (math.abs(el[5] - vk(h)) > param.eps_v) {
        if (t + dt > t1 && t < t1) {
            dt = t1 - t;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i)
            dt = param.step - dt;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i);
            dt = param.step
        }
        if (t + dt > t2 && t < t2) {
            dt = t2 - t;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i)
            dt = param.step - dt;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i)
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
        t += dt;
        if (el[5] < vk(h) && dt === param.step) {
            html(t, el, i);
        }
        if (math.abs(el[5] - vk(h)) < param.eps_v) {
            html(t, el, i);
        }
        if (el[5] > vk(h)) {
            el = el_prev.slice();
            t = t_prev;
            el_prev = el.slice();
            t_prev = t;
            dt = dt / 10;
        }
    }
    i++
}

export {
    html,
    render,
    consoleLog
}